import React, { createContext, useContext, useState, useEffect, useCallback } from 'react';

export interface SessionData {
  id: string;
  name: string;
  timestamp: number;
  vcfPath?: string;
  selectedRegion?: string;
  selectedGenotypes?: string[];
  publicVcfFolder?: string;
  cohortFolder?: string;
  mode?: 'individual' | 'cohort';
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

  // Save sessions to localStorage whenever they change
  useEffect(() => {
    try {
      localStorage.setItem(STORAGE_KEY, JSON.stringify(sessions));
    } catch (error) {
      console.error('Error saving sessions to localStorage:', error);
    }
  }, [sessions]);

  // Save current session to localStorage whenever it changes
  useEffect(() => {
    try {
      if (currentSession) {
        localStorage.setItem(CURRENT_SESSION_KEY, JSON.stringify(currentSession));
      } else {
        localStorage.removeItem(CURRENT_SESSION_KEY);
      }
    } catch (error) {
      console.error('Error saving current session to localStorage:', error);
    }
  }, [currentSession]);

  const saveSession = useCallback((session: Omit<SessionData, 'id' | 'timestamp'>): string => {
    const newSession: SessionData = {
      ...session,
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
      return { ...prev, ...updates };
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




