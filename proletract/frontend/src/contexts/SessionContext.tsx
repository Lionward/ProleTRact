import React, { createContext, useContext, useState, useEffect, useCallback } from 'react';

export type SessionMode = 'individual' | 'cohort-read' | 'cohort-assembly';

export interface SessionData {
  id: string;
  name: string;
  timestamp: number;
  vcfPath?: string;
  selectedRegion?: string;
  selectedGenotypes?: string[];
  publicVcfFolder?: string;
  cohortFolder?: string;
  mode?: SessionMode;
  notes?: string;
}

interface SessionContextType {
  sessions: SessionData[];
  currentSession: SessionData | null;
  saveSession: (session: Omit<SessionData, 'id' | 'timestamp'>) => string;
  loadSession: (sessionId: string) => void;
  deleteSession: (sessionId: string) => void;
  updateCurrentSession: (updates: Partial<SessionData>) => void;
  clearCurrentSession: () => void;
}

const SessionContext = createContext<SessionContextType | undefined>(undefined);

const STORAGE_KEY = 'proletract_sessions';
const CURRENT_SESSION_KEY = 'proletract_current_session';

/** Extract only serializable fields to avoid "Converting circular structure to JSON" when DOM/React refs leak in. */
function sanitizeSessionData<T extends Partial<SessionData>>(raw: T): T {
  return {
    ...(raw.id !== undefined && { id: String(raw.id) }),
    ...(raw.name !== undefined && { name: String(raw.name) }),
    ...(raw.timestamp !== undefined && { timestamp: Number(raw.timestamp) }),
    ...(raw.vcfPath !== undefined && { vcfPath: String(raw.vcfPath ?? '') }),
    ...(raw.selectedRegion !== undefined && { selectedRegion: String(raw.selectedRegion ?? '') }),
    ...(raw.selectedGenotypes !== undefined && {
      selectedGenotypes: Array.isArray(raw.selectedGenotypes)
        ? raw.selectedGenotypes.map((g) => String(g)).filter(Boolean)
        : [],
    }),
    ...(raw.publicVcfFolder !== undefined && { publicVcfFolder: String(raw.publicVcfFolder ?? '') }),
    ...(raw.cohortFolder !== undefined && { cohortFolder: String(raw.cohortFolder ?? '') }),
    ...(raw.mode !== undefined && {
      mode: raw.mode === 'cohort-read' || raw.mode === 'cohort-assembly' ? raw.mode : ('individual' as SessionMode),
    }),
    ...(raw.notes !== undefined && { notes: raw.notes != null ? String(raw.notes) : undefined }),
  } as T;
}

/** Safely stringify to avoid "Converting circular structure to JSON" from DOM/React refs. */
function safeStringify(obj: unknown): string {
  const seen = new WeakSet();
  return JSON.stringify(obj, (_, value) => {
    if (value != null && typeof value === 'object') {
      if (seen.has(value)) return undefined;
      seen.add(value);
    }
    if (typeof value === 'function' || (typeof HTMLElement !== 'undefined' && value instanceof HTMLElement)) {
      return undefined;
    }
    return value;
  });
}

export const SessionProvider: React.FC<{ children: React.ReactNode }> = ({ children }) => {
  const [sessions, setSessions] = useState<SessionData[]>([]);
  const [currentSession, setCurrentSession] = useState<SessionData | null>(null);

  // Load sessions from localStorage on mount
  useEffect(() => {
    try {
      const storedSessions = localStorage.getItem(STORAGE_KEY);
      if (storedSessions) {
        setSessions(JSON.parse(storedSessions));
      }

      const storedCurrent = localStorage.getItem(CURRENT_SESSION_KEY);
      if (storedCurrent) {
        setCurrentSession(JSON.parse(storedCurrent));
      }
    } catch (error) {
      console.error('Error loading sessions from localStorage:', error);
    }
  }, []);

  // Save sessions to localStorage whenever they change (sanitize + safeStringify to avoid circular refs)
  useEffect(() => {
    try {
      const safe = sessions.map(s => sanitizeSessionData(s) as SessionData);
      localStorage.setItem(STORAGE_KEY, safeStringify(safe));
    } catch (error) {
      console.error('Error saving sessions to localStorage:', error);
    }
  }, [sessions]);

  // Save current session to localStorage whenever it changes (sanitize + safeStringify to avoid circular refs)
  useEffect(() => {
    try {
      if (currentSession) {
        const safe = sanitizeSessionData(currentSession) as SessionData;
        localStorage.setItem(CURRENT_SESSION_KEY, safeStringify(safe));
      } else {
        localStorage.removeItem(CURRENT_SESSION_KEY);
      }
    } catch (error) {
      console.error('Error saving current session to localStorage:', error);
    }
  }, [currentSession]);

  const saveSession = useCallback((session: Omit<SessionData, 'id' | 'timestamp'>): string => {
    const sanitized = sanitizeSessionData(session) as Omit<SessionData, 'id' | 'timestamp'>;
    const newSession: SessionData = {
      ...sanitized,
      id: `session_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`,
      timestamp: Date.now(),
    };

    setSessions(prev => {
      // Remove old session with same name if it exists
      const filtered = prev.filter(s => s.name !== newSession.name);
      return [...filtered, newSession].sort((a, b) => b.timestamp - a.timestamp);
    });

    setCurrentSession(newSession);
    return newSession.id;
  }, []);

  const loadSession = useCallback((sessionId: string) => {
    const session = sessions.find(s => s.id === sessionId);
    if (session) {
      setCurrentSession(session);
    }
  }, [sessions]);

  const deleteSession = useCallback((sessionId: string) => {
    setSessions(prev => prev.filter(s => s.id !== sessionId));
    if (currentSession?.id === sessionId) {
      setCurrentSession(null);
    }
  }, [currentSession]);

  const updateCurrentSession = useCallback((updates: Partial<SessionData>) => {
    setCurrentSession(prev => {
      if (!prev) return null;
      const safePrev = sanitizeSessionData(prev) as SessionData;
      const safeUpdates = sanitizeSessionData(updates);
      return { ...safePrev, ...safeUpdates } as SessionData;
    });
  }, []);

  const clearCurrentSession = useCallback(() => {
    setCurrentSession(null);
  }, []);

  return (
    <SessionContext.Provider
      value={{
        sessions,
        currentSession,
        saveSession,
        loadSession,
        deleteSession,
        updateCurrentSession,
        clearCurrentSession,
      }}
    >
      {children}
    </SessionContext.Provider>
  );
};

export const useSession = (): SessionContextType => {
  const context = useContext(SessionContext);
  if (!context) {
    throw new Error('useSession must be used within a SessionProvider');
  }
  return context;
};




